import np_config
import pathlib

zk_config = np_config.from_zk("/projects/np_datajoint/defaults/configuration")
datajoint_exp = pathlib.Path(zk_config["sorting"]["local_inbox"]) / 'ks_paramset_idx_1'
metrics_path = pathlib.Path("1187668018_614608_20220628/1187668018_614608_20220628_probeA_sorted/continuous/Neuropix-PXI-100.0/metrics.csv")
np_exp = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')
fulldj_path = (datajoint_exp / metrics_path).resolve().absolute()
fulldj_path.relative_to(np_exp)