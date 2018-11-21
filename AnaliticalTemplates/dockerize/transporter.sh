#!/usr/bin/env bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-c|--container)
	    CONTAINER="$2"
	    shift
	    shift
	    ;;
	-s|--source)
	    SPATH="$2"
	    shift
	    shift
	    ;;
	-p|--path)
	    OUTPATH="$2"
	    shift
	    shift
	    ;;
	-b|--builder)
            BUILDER="$2"
            shift
            shift
            ;;
    esac
done

set -- "${POSITIONAL[@]}"

EPOINT="${BUILDER} ${SPATH} ${OUTPATH}"

if [[ "$CONTAINER" = "centos" ]]; then
    CMD="docker run --rm -it -v /Users/enrico/:/userhome --user $(id -u) rootframework:centos ${EPOINT}"
    ${CMD}
elif [[ "$CONTAINER" = "ubuntu" ]]; then
    CMD="docker run --rm -it -v /Users/enrico/:/userhome --user $(id -u) rootframework:ubuntu ${EPOINT}"
    ${CMD}
else
    echo "Unkown container selected"
fi

