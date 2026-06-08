#!/home/dell/miniconda3/envs/PathoSource/bin/python
import multiprocessing
import os
import sys

from pathosource_refactor.config import build_argument_parser, build_pipeline_config
from pathosource_refactor.context import RuntimeContext, set_runtime_context
from pathosource_refactor.workflow import main_process

__author__ = "wsh"
__version__ = "2.0.0"
__date__ = "20260330"


def bootstrap_runtime():
    parser = build_argument_parser()
    argv = parser.parse_args()
    command_line = " ".join(sys.argv)
    print(command_line)
    sys.stdout.flush()

    pipeline_config = build_pipeline_config(argv, multiprocessing.cpu_count())

    print("开始数据分析")
    print(pipeline_config.config_out)
    sys.stdout.flush()

    os.makedirs(pipeline_config.output_dir, exist_ok=True)
    os.chdir(pipeline_config.output_dir)

    set_runtime_context(
        RuntimeContext(
            input_path=pipeline_config.input_path,
            run_path=pipeline_config.output_dir,
            species=pipeline_config.species,
            threads=pipeline_config.threads,
            meta=pipeline_config.meta,
            ref=pipeline_config.ref,
            mode=pipeline_config.mode,
            cgmlstana=pipeline_config.cgmlstana,
            gubbins=pipeline_config.gubbins,
            msamethod=pipeline_config.msamethod,
            treemethod=pipeline_config.treemethod,
            bootstrap=pipeline_config.bootstrap,
            mltype=pipeline_config.mltype,
            cgmlstversion=pipeline_config.cgmlstversion,
            speciesdb=pipeline_config.speciesdb,
            fametadb=pipeline_config.fametadb,
            resources=pipeline_config.resources,
        )
    )
    return pipeline_config


def run() -> None:
    pipeline_config = bootstrap_runtime()
    main_process(pipeline_config.meta, pipeline_config.ref, pipeline_config.input_path)


if __name__ == "__main__":
    run()
