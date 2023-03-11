from pathlib import Path
import pandas as pd
import os

def create_path() -> None:
    base_path = Path('//allen/programs/mindscope/workgroups/np-behavior')
    full_path = Path(base_path, 'tissuecyte', '615047', 'probe_annotations_615047.csv')
    return full_path

def exists(path: Path) -> None:
    print(path.exists())

def parent() -> None:
    metrics_path = Path('//allen/programs/mindscope/workgroups/np-exp/1178173272_608671_20220518/1178173272_608671_20220518_probeB_sorted/continuous/Neuropix-PXI-100.0/metrics.csv')
    print(metrics_path.parent)
    spike_times = Path(metrics_path.parent, 'spike_times.npy')
    print(spike_times.exists())

def iter() -> None:
    p = Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte')
    print([x for x in p.iterdir() if x.is_dir()])


def join() -> None:
    base_path = Path('//allen/programs/mindscope/workgroups/np-behavior')
    full_path = base_path / 'tissuecyte' / '608671'
    join_path = base_path.joinpath('tissuecyte', '608671')
    print(type(full_path), full_path)
    print(full_path.exists())
    print(join_path.exists())

def glob_lib(path: Path) -> None:
    if path.is_dir():
        for p in path.glob('Probe_F2*'):
            print(p.name)

def rglob_lib(path: Path) -> None:
    if path.is_dir():
        for p in path.rglob('Probe_F2*'):
            print(p.name)

if __name__ == '__main__':
    str_path = r'/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/615047/probe_annotations_615047.csv'
    #os_path = os.path.join('/allen/programs/mindscope/workgroups/np-behavior', 'tissuecyte', '615047', 'probe_annotations_615047.csv')

    #pathlib_path = Path('//allen/programs/mindscope/workgroups/np-behavior/tissuecyte/615047/probe_annotations_615047.csv')
    pathlib_path = create_path()

    annotation_file_str = pd.read_csv(str_path, index_col=0)
    #os_file = pd.read_csv(os_path, index_col=0)
    annotation_file_pathlib = pd.read_csv(pathlib_path, index_col=0)
    print('Str', annotation_file_str)
    print('Pathlib', annotation_file_pathlib)
