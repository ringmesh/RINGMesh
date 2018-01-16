travis_job_id=$1
travis_job_number=$2

echo "Test"
echo "$travis_job_id"
echo "$travis_job_number"
  
curl -X POST \
  --data-urlencode "payload={ \
  \"channel\": \"#travis\", \
  \"username\": \"Travis CI\", \
  \"text\": \"Problem on master: see <https://travis-ci.org/ringmesh/RINGMesh/jobs/329406106|Failed job> $travis_job_id and $travis_job_number.\"}" \
  https://hooks.slack.com/services/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
