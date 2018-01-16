body='{
"request": {
  "branch":"master"
}}'

  
curl -X POST \
  --data-urlencode "payload={ \
  \"channel\": \"#travis\", \
  \"username\": \"Travis CI\", \
  \"on_success\": \"never\", \
  \"on_failure\": \"always\", \
  \"text\": \"This is posted to #travis channel and comes from a bot named Travis CI.\", \
  \"icon_emoji\": \":x:\"}" \
  https://hooks.slack.com/services/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
