"""Functions for running the inDrop pipeline on the cloud."""

import json
import subprocess
import re
import os
import copy
import tempfile
import time
import yaml

from jinja2 import Environment, PackageLoader, select_autoescape

from genometools import gcloud

from .. import __version__, _root
from . import config

_GSURL_PAT = re.compile(r'gs://(.*?)/(.*)')

_TEMPLATE_ENV = Environment(
    loader=PackageLoader('singlecell',
                         os.path.join('data', 'templates')),
    autoescape=select_autoescape(['html', 'xml']))


def run_subprocess(cmd):
    subproc = subprocess.Popen(cmd, shell=True)
    subproc.communicate()
    return subproc.returncode


def split_gsurl(gsurl):
    import re
    match = _GSURL_PAT.match(gsurl)
    return match.group(1), match.group(2)


def upload_singlecell_package(conf, overwrite=True):

    gc = conf['gcloud']
    
    # create source tarball
    package_dir = os.path.abspath(os.path.join(_root, '..'))
    cmd = 'cd "%s" && python setup.py sdist' % package_dir
    run_subprocess(cmd)
    
    # upload it to the cloud
    with open(gc['key_file']) as fh:
        key = json.load(fh)

    project_id = gc['project_id']
    client = gcloud.storage.get_client(key, project_id)
    package_file = os.path.join(package_dir, 'dist', 'singlecell-%s.tar.gz'
                                % __version__)
    assert os.path.isfile(package_file)
    bucket, path = split_gsurl(gc['singlecell_package_file'])
    gcloud.storage.upload_file(client, bucket, package_file, path,
                               overwrite=overwrite)


def run_pipeline(config_file):

    conf, errors = config.read_config(config_file)

    input_ = conf['input']
    output = conf['output']
    params = conf['parameters']
    pipeline = conf['pipeline']
    gc = conf['gcloud']
    
    template = _TEMPLATE_ENV.get_template(os.path.join('gcloud', 'indrop_pipeline_wrapper.sh'))
    
    timestamp = time.strftime('%Y-%m-%d-%H-%M-%S')
    
    if pipeline['skip_read_processing'] and not pipeline['skip_mapping']:
        raise ValueError('If you skip processing, you must also skip mapping!')
        
    # create copy of config file for the Google Compute instance
    conf2 = copy.deepcopy(conf)
    conf2['input']['barcode_read_file'] = '/root/data/barcode_reads.fastq.gz'
    conf2['input']['mrna_read_file'] = '/root/data/mrna_reads.fastq.gz'
    conf2['input']['genome_file'] = '/root/data/genome.fa.gz'
    conf2['input']['genome_annotation_file'] = '/root/data/genome_annotations.gtf.gz'
    conf2['input']['star_index_dir'] = '/root/data/star_index'
    conf2['output']['output_dir'] = '/root/output'
    conf2['output']['allow_nonempty_output_dir'] = True
    conf2['parameters']['use_docker'] = True

    # upload config file to the cloud
    with open(gc['key_file']) as fh:
        key = json.load(fh)    

    config_file_name = 'config_%s.yaml' % timestamp
    config_gsurl = gc['pipeline_config_dir'] + '/' + config_file_name
    with tempfile.NamedTemporaryFile(mode='w+', encoding='utf-8') as tfh:
        # write modified configuration file to disk
        tfh.write(yaml.dump(conf2, indent=4, default_flow_style=False))
        tfh.flush()
        
        # upload the file to Google Cloud Storage
        project_id = gc['project_id']
        client = gcloud.storage.get_client(key, project_id)
        
        bucket, path = split_gsurl(config_gsurl)
        gcloud.storage.upload_file(client, bucket, tfh.name, path,
                                   overwrite=False)
    
    startup_script = template.render(
        script_dir=gc['script_dir'],
        barcode_read_file=input_['barcode_read_file'],
        mrna_read_file=input_['mrna_read_file'],
        genome_file=input_['genome_file'],
        genome_annotation_file=input_['genome_annotation_file'],
        index_dir=input_['star_index_dir'],
        output_dir=output['output_dir'],
        singlecell_package_file=gc['singlecell_package_file'],
        singlecell_package_version=__version__,
        pipeline_config_file=config_gsurl,
        skip_read_processing=pipeline['skip_read_processing'],
        skip_mapping=pipeline['skip_mapping'],
        skip_aligned_read_processing=pipeline['skip_aligned_read_processing'],
        skip_expression_quantification=pipeline['skip_expression_quantification'],
        skip_qc_plot_generation=pipeline['skip_qc_plot_generation'],
        self_destruct=True,
    )
    
    project_id = gc['project_id']
    zone = gc['zone']
    machine_type = gc['machine_type']
    disk_size_gb = gc['disk_size_gb']
    instance_config = gcloud.compute.InstanceConfig(
        project_id, zone,
        machine_type=machine_type,
        disk_size_gb=disk_size_gb)

    instance_name = 'indrop-pipeline-%s' % timestamp
    
    with open(gc['key_file']) as fh:
        key = json.load(fh)
    cred = gcloud.get_compute_credentials(key)
    
    op_name = instance_config.create_instance(
        cred, instance_name, startup_script=startup_script,
        wait_until_done=True)

    return instance_name

"""
def run_pipeline(config_file):

    conf = config.read_config(config_file)

    input_ = conf['input']
    params = conf['parameters']
    output_dir = conf['output_dir']
    
    script_dir = conf['gcloud']['script_dir']
    script_dir = script_dir.rstrip('/')

    with open(conf['gcloud']['key_file']) as fh:
        key = json.load(fh)

    project_id = key['project_id']

    #cred = gcloud.get_compute_credentials(key)
    
    client = gcloud.storage.get_client(key, project_id)
    gcloud.tasks.upload_scripts(client, script_dir)
"""