#!/bin/bash


docker build . -t masq-development:latest

docker run --rm -v $(pwd):/wd -it masq-development:latest \
    /bin/bash -c "cd /wd; ./integration_test.sh"
