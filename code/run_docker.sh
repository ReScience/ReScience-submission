#!/bin/bash
set -eo pipefail

docker build -t python_tabak -f Dockerfile .
docker run -i -v $(pwd):/home/docker/code -v $(pwd)/../article:/home/docker/article -t python_tabak
