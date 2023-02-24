# generates a csv with the paths for the metrics file

import pandas as pd
import argparse
import os
import pathlib

parser = argparse.ArgumentParser()
parser.add_argument('--mouseID', help='Mouse ID of session', required=True)

# gets the directory where the metrics file is
def get_metrics_directory(base_path: str,  mouse_id: str):
    directories = os.listdir(base_path)
    probe_directories = []

    for d in directories:
        if mouse_id in d:
            probe_directories.append(d)
    
    return probe_directories

def generate_templeton_metric_path_days(mouse_id: str, base_path: str, record_node: str):
    record_node = 'Record Node {}'.format(record_node)
    base_path = pathlib.Path(base_path)
    mouse_dirs = get_metrics_directory(base_path, mouse_id)
    probe_metrics_dirs = {}
    metrics_path_days = {}

    for directory in mouse_dirs:
        date = directory[0:directory.index('_')]
        probe_dirs = [d for d in os.listdir(os.path.join(base_path, directory, record_node, 'experiment1', 'recording1', 'continuous')) 
                      if os.path.isdir(os.path.join(base_path, directory, record_node, 'experiment1', 'recording1', 'continuous'))]
        probe_metrics_dirs [date] = [os.path.join(base_path, directory, record_node, 'experiment1', 'recording1', 'continuous', d) for d in probe_dirs]
            
    for date in probe_metrics_dirs:
        for d in probe_metrics_dirs[date]:
            files = [os.path.join(d, f) for f in os.listdir(d)]
            for f in files:
                if 'waveform_metrics.csv' in f:
                    waveform_metrics_path = f
                    if date not in metrics_path_days:
                        metrics_path_days[date] = [waveform_metrics_path]
                    else:
                        metrics_path_days[date].append(waveform_metrics_path)

    return metrics_path_days

# gets the path for the metrics csv
def generate_metrics_path_days(base_path, mouse_id):
    mouse_dirs = get_metrics_directory(base_path, mouse_id)
    probe_metrics_dirs = []

    for directory in mouse_dirs:
        probe_dirs = [d for d in os.listdir(os.path.join(base_path, directory)) if os.path.isdir(os.path.join(base_path, directory, d))]

        probe_metrics_dirs.append([os.path.join(base_path, directory, d, 'continuous') for d in probe_dirs if os.path.exists(os.path.join(base_path, directory, d, 'continuous'))])
    
    metrics_path_days = {}
    
    for directory in probe_metrics_dirs:
        if len(directory) > 0:
            for d in directory:
                date = d[d.index('_')+1:]
                date = date[date.index('_')+1:date.index('\\')]

                files = [os.path.join(d, f) for f in os.listdir(d)]
    
                for f in files:
                    metrics = os.listdir(os.path.join(d, f))

                    waveform_metrics_path = [os.path.join(d, f, m) for m in metrics if m == 'metrics.csv']

                    if len(waveform_metrics_path) > 0:
                        waveform_metrics_path = waveform_metrics_path[0]
                        if date not in metrics_path_days:
                            metrics_path_days[date] = [waveform_metrics_path]
                        else:
                            metrics_path_days[date].append(waveform_metrics_path)

    #print(metrics_path_days)
    return metrics_path_days

if __name__ == '__main__':
    base_path = pathlib.Path('//allen/programs/mindscope/workgroups/np-exp')
    output_path = pathlib.Path('//allen/programs/mindscope/workgroups/np-behavior')

    args = parser.parse_args()
    mouse_id = args.mouseID

    #generate_metrics_path_days(base_path, output_path, mouse_id)
    generate_templeton_metric_path_days(mouse_id)