#!/bin/bash

travis_job_id=$1
travis_branch=$2
  
curl -X POST \
  --data-urlencode "payload={ \
  \"channel\": \"#travis\", \
  \"username\": \"Travis CI\", \
  \"text\": \"Problem on $travis_branch branch: see <https://travis-ci.org/ringmesh/RINGMesh/jobs/$travis_job_id|Failed job on Travis>.\"}" \
  https://hooks.slack.com/services/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
