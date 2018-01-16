travis_job_id=$1
travis_job_number=$2

echo "Test"
echo "$travis_job_id"
echo "$travis_job_number"
  
curl -X POST \
  --data-urlencode "payload={ \
  \"channel\": \"#travis\", \
  \"username\": \"Travis CI\", \
  \"on_success\": \"change\", \
  \"on_failure\": \"always\", \
  \"text\": \"Problem with $travis_job_id and $travis_job_number.\", \
  \"icon_emoji\": \":x:\"}" \
  https://hooks.slack.com/services/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
