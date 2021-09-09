TIMESTEP=300.
OVERWRITE=false


init:
	@declare -a opts=();\
	opts+=('forecast');\
	opts+=('../hgrid.ll');\
	opts+=(--timestep=${TIMESTEP});\
	if [[ ${OVERWRITE} = true ]];\
	then \
		opts+=("--overwrite");\
	fi;\
	opts+=(--forecast-days=2);\
	opts+=(--spinup-days=1);\
	opts+=(--major-constituents);\
	opts+=(--iettype-3);\
	opts+=(--ifltype-3);\
	opts+=(--sflux gfs);\
	opts+=(--modules-init=/usr/local/Modules/default/init/sh);\
	opts+=(--modulepath="/sciclone/home20/jrcalzada/.local/Modules/modulefiles:/usr/local/Modules/modulefiles:/usr/local/Modules/default/modulefiles");\
	opts+=(--module slurm);\
	opts+=(--module schism);\
	opts+=(--log-level=info);\
	opts+=(slurm);\
	opts+=(--walltime=10m);\
	opts+=(--ntasks=30);\
	echo "pyschism forecast init" "$${opts[@]}";\
	/sciclone/home20/jrcalzada/.miniconda3/envs/pyschism/bin/pyschism forecast init "$${opts[@]}"
