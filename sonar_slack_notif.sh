#!/bin/bash
  
echo Script_sonar_OK  

curl -X POST \
  --data-urlencode "payload={ \
  \"channel\": \"#general\", \
  \"username\": \"Sonar\", \
  \"text\": \"Coverage: see https://sonarcloud.io/api/badges/measure?key=ringmesh&metric=coverage and https://sonarcloud.io/dashboard/index/ringmesh\"}" \
https://hooks.slack.com/services/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
