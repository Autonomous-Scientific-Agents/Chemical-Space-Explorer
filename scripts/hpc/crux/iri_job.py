"""Write an IRI JobSpec (ALCF IRI API v1) for Crux; submit it with the Genesis iri-api skill.

    python scripts/hpc/crux/iri_job.py setup --user USER --project PROJ > setup.json
    python scripts/hpc/crux/iri_job.py run --user USER --project PROJ --name cnof12 --nodes 2 \
        -- --kind exact --max-atoms 12 --elements C N O F --shards 4000 > run.json

    # genesis-skills/skills/amsc-skills/iri-api:
    python3 scripts/iri_api_call.py --facility alcf call --operation-id launchJob \
        --path-param resource_id=8b9b42f7-572a-4909-8472-a0453436304c --json-file run.json --ensure-token

Paths are absolute: IRI passes $HOME and $USER literally in JSON fields.
Outputs: small files in /home/USER/chemuniverse-runs/NAME (readable with the IRI `view`
operation); large shard data in /lus/grand/projects/PROJ/chemuniverse/NAME.
"""
import argparse
import json
import shlex
import sys
import time

CRUX = "8b9b42f7-572a-4909-8472-a0453436304c"
REPO = "https://github.com/Autonomous-Scientific-Agents/Chemical-Universe.git"
PROXY = "http://proxy.alcf.anl.gov:3128"


def spec(a, name, command):
    home = f"/home/{a.user}"
    out = f"{home}/chemuniverse-runs/{name}"
    return {
        "name": f"chemu-{name}"[:15],
        "executable": "/bin/bash",
        "arguments": ["-c", f"mkdir -p {out} && {command}"],
        "directory": home,
        "stdout_path": f"{out}/job.out",
        "stderr_path": f"{out}/job.err",
        "resources": {"node_count": a.nodes},
        "attributes": {
            "duration": a.walltime * 60,
            "queue_name": a.queue,
            "account": a.project,
            "custom_attributes": {"filesystems": "home:grand"},
        },
    }


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("kind", choices=["setup", "run", "probe"],
                   help="setup: build the environment; run: a sharded job; probe: report how the job is launched")
    p.add_argument("--user", required=True, help="ALCF username")
    p.add_argument("--project", required=True, help="ALCF project with a Crux allocation")
    p.add_argument("--name", default=None, help="run name (default: setup-<timestamp>)")
    p.add_argument("--nodes", type=int, default=1)
    p.add_argument("--queue", default="debug", help="debug (1-4 nodes, 2 h), workq-route (<=184 nodes, 24 h), "
                                                     "preemptable (<=10 nodes, 72 h)")
    p.add_argument("--walltime", type=int, default=60, help="minutes")
    p.add_argument("--prefix", default=None, help="environment prefix (default /home/USER/chemuniverse-crux)")
    p.add_argument("--run-dir", default=None, help="default /lus/grand/projects/PROJECT/chemuniverse/NAME")
    argv = sys.argv[1:]
    plan_args = argv[argv.index("--") + 1:] if "--" in argv else []
    a = p.parse_args(argv[:argv.index("--")] if "--" in argv else argv)
    prefix = a.prefix or f"/home/{a.user}/chemuniverse-crux"
    env = f"export CHEMUNIVERSE_PREFIX={prefix} https_proxy={PROXY} http_proxy={PROXY}"
    if a.kind == "probe":
        name = a.name or time.strftime("probe-%Y%m%d-%H%M%S")
        command = ("echo host=$(hostname) rank=${PALS_RANKID:-${PMI_RANK:-none}} nproc=$(nproc); "
                   "echo nodefile:; cat $PBS_NODEFILE; grep -m1 'model name' /proc/cpuinfo; "
                   "which mpiexec python3; python3 --version; ls -d /lus/grand/projects/" + a.project)
    elif a.kind == "setup":
        name = a.name or time.strftime("setup-%Y%m%d-%H%M%S")
        src = f"{prefix}/src"
        command = (f"{env} && mkdir -p {prefix} && "
                   f"(test -d {src}/.git || git clone -q {REPO} {src}) && git -C {src} pull -q --ff-only && "
                   f"bash {src}/scripts/hpc/crux/setup_env.sh")
    else:
        if not a.name:
            p.error("run needs --name")
        if not plan_args:
            p.error("run needs shard plan arguments after --, e.g. -- --kind exact --max-atoms 12 ...")
        name = a.name
        run_dir = a.run_dir or f"/lus/grand/projects/{a.project}/chemuniverse/{name}"
        summary = f"/home/{a.user}/chemuniverse-runs/{name}"
        command = (f"{env} && bash {prefix}/src/scripts/hpc/crux/run_shards.sh {shlex.quote(run_dir)} "
                   f"{shlex.quote(summary)} -- {' '.join(shlex.quote(x) for x in plan_args)}")
    json.dump(spec(a, name, command), sys.stdout, indent=2)
    print()


if __name__ == "__main__":
    main()
