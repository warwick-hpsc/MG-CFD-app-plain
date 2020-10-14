#!/bin/bash

## There are many macro flags that control compilation paths.
## To ensure all paths compile, test all permutations:

_blank="-"

_flux_options=("$_blank" "FLUX_CRIPPLE" "FLUX_PRECOMPUTE_EDGE_WEIGHTS")
_simd=("$_blank" "SIMD")
_simd_options=("$_blank" "COLOUR" "MANUAL" "MANUAL_GATHER" "MANUAL_SCATTER")
_misc_options=("$_blank" "TIME" "PAPI" "OMP" "OMP_SCATTERS")

_compile() (
	set -e
	set -u

	for _f in ${_flux_options[@]}; do
		if [ "$_f" = "-" ]; then
			bf1=""
		else
			bf1="-D$_f "
		fi

		for _s in ${_simd[@]}; do
			bf2="$bf1"
			if [ "$_s" != "$_blank" ]; then
				bf2+="-D${_s} -DDBLS_PER_SIMD=2 "
			fi

			for _so in ${_simd_options[@]}; do
				bf3="$bf2"
				if [ "$_so" != "$_blank" ]; then
					bf3+="-D${_so} "
				fi

				if [ "$_so" = "COLOUR" ]; then
					_colour_options=("$_blank" "BIN_COLOURED_VECTORS" "BIN_COLOURED_CONTIGUOUS")
				else
					_colour_options=("$_blank")
				fi
		
				for _c in ${_colour_options[@]}; do 
					bf4="$bf3"
					if [ "$_c" != "$_blank" ]; then
						bf4+="-D${_c} "
					fi

					for _m in ${_misc_options[@]}; do
						bf5="$bf4"
						if [ "$_m" != "$_blank" ]; then
							bf5+="-D${_m} "
						fi

						bf_final="$bf5"
						echo "Testing build flags: $bf_final"
						export BUILD_FLAGS="$bf_final"
						bf4b=`echo "$bf_final" | sed -s "s/ //g"`
						bf4b+="-DINSN_SET=Host"
						bin_filename="euler3d_cpu_double_${COMPILER}${bf4b}.b"
						make -j`nproc`
						if [ ! -f bin/"${bin_filename}" ]; then
							## Run make again to print error clearly in stdout:
							make
							echo "ERROR: compilation failed to generate bin: ${bin_filename}"
							echo "       $bf_final"
							exit 1
						fi
					done
				done
			done
		done
	done
)

_validate_binaries() (
	set -e

	cd ~/Datasets/Live/M6_wing/rodinia.regenerated
	ls ~/Working.Github.staging/MG-CFD-app-plain.SIMD/bin | while read F ; do
		echo "Checking: $F"
		~/Working.Github.staging/MG-CFD-app-plain.SIMD/bin/"$F" -i input.dat -g 5 -v -p papi.conf
	done
)

_compile

_validate_binaries
