"""Download ISIMIP3b GGCMI yields (rainfed, future SSPs), aggregate to country annual means
weighted by 2015 harvested area (ISIMIP 2015soc landuse-15crops) x fractional country mask."""
import json, os, re, sys, time, tempfile, requests, threading
import numpy as np, pandas as pd, xarray as xr
from concurrent.futures import ThreadPoolExecutor, as_completed
D = sys.argv[1]  # prep dir
OUT = os.path.join(D, "out"); os.makedirs(OUT, exist_ok=True)
ISO = ["BFA","BDI","COM","CIV","COD","ETH","GHA","KEN","MWI","MLI","NER","NGA","RWA","SEN","SLE","ZAF","UGA","TZA","ZMB"]
LAT = slice(25, -36); LON = slice(-20, 53)
LU = {"mai":"maize","ri1":"rice","ri2":"rice","soy":"oil_crops_soybean","swh":"temperate_cereals","wwh":"temperate_cereals",
      "mil":"tropical_cereals","sor":"tropical_cereals","bea":"pulses","pot":"temperate_roots","cas":"tropical_roots"}
masks = xr.open_dataset(os.path.join(D,"countrymasks-fractional_30arcmin.nc")).sel(lat=LAT, lon=LON)
lu = xr.open_dataset(os.path.join(D,"landuse-15crops_2015soc_annual_2015_2100.nc"), decode_times=False).isel(time=0).sel(lat=LAT, lon=LON)
W = {}
for c, v in LU.items():
    area = (lu[f"{v}_rainfed"].fillna(0) + lu[f"{v}_irrigated"].fillna(0)).values
    W[c] = {iso: area * masks[f"m_{iso}"].fillna(0).values for iso in ISO}
recs = json.load(open(os.path.join(D, "yield_datasets.json")))
sel = [r for r in recs if r["scen"] in ("ssp126","ssp585") and r["irr"]=="noirr" and r["crop"] in ("mai","ri1","soy","swh","mil","sor","bea","pot","cas")]
done_f = os.path.join(OUT, "done.txt"); done = set(open(done_f).read().split()) if os.path.exists(done_f) else set()
sel = [r for r in sel if r["name"] not in done]
print("todo", len(sel), flush=True)
lock = threading.Lock()
def work(r):
    path, size, url = r["files"][0]
    y0, y1 = map(int, re.search(r"_(\d{4})_(\d{4})\.nc$", path).groups())
    fd, tmp = tempfile.mkstemp(suffix=".nc", dir=OUT); os.close(fd)
    for a in range(4):
        try:
            with requests.get(url, stream=True, timeout=120) as resp:
                resp.raise_for_status()
                with open(tmp, "wb") as f:
                    for ch in resp.iter_content(1 << 20): f.write(ch)
            break
        except Exception as e:
            if a == 3: raise
            time.sleep(5 * (a + 1))
    try:
        ds = xr.open_dataset(tmp, decode_times=False)
        var = [v for v in ds.data_vars if v.startswith("yield")][0]
        y = ds[var].sel(lat=LAT, lon=LON).values.astype("float64")  # t/ha dry matter
        ds.close()
    finally:
        os.remove(tmp)
    y[(y > 1e19) | (y < 0)] = np.nan
    years = np.arange(y0, y1 + 1)[: y.shape[0]]
    rows = []
    for iso in ISO:
        w = W[r["crop"]][iso]
        valid = np.isfinite(y)
        num = np.nansum(np.where(valid, y * w, 0), axis=(1, 2)); den = np.sum(np.where(valid, w, 0), axis=(1, 2))
        # fallback: unweighted mean over country cells when crop area is zero
        m = masks[f"m_{iso}"].fillna(0).values
        num2 = np.nansum(np.where(valid, y * m, 0), axis=(1, 2)); den2 = np.sum(np.where(valid, m, 0), axis=(1, 2))
        for i, yr in enumerate(years):
            rows.append((r["model"], r["gcm"], r["scen"], r["sens"], r["crop"], iso, int(yr),
                         num[i] / den[i] if den[i] > 0 else np.nan, num2[i] / den2[i] if den2[i] > 0 else np.nan, den[i]))
    df = pd.DataFrame(rows, columns=["model","gcm","scen","co2","crop","iso3","year","yield_areaw","yield_unw","crop_area_ha_weight"])
    with lock:
        f = os.path.join(OUT, "isimip3b_country_annual.csv")
        df.to_csv(f, mode="a", header=not os.path.exists(f), index=False)
        open(done_f, "a").write(r["name"] + "\n")
    return r["name"]
n = 0; t0 = time.time()
with ThreadPoolExecutor(8) as ex:
    futs = [ex.submit(work, r) for r in sel]
    for fu in as_completed(futs):
        n += 1
        try: fu.result()
        except Exception as e: print("FAIL", e, flush=True)
        if n % 20 == 0: print(n, len(sel), f"{time.time()-t0:.0f}s", flush=True)
print("DONE", flush=True)
