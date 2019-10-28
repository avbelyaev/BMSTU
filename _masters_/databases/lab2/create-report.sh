#!/usr/bin/env bash

REPORT_NAME=report-2-belyaev.md

echo "removing old report"
rm $REPORT_NAME

echo "creating report ${REPORT_NAME}"
touch $REPORT_NAME

# copy md
cat q1-2.md >> $REPORT_NAME

# copy SQLs
echo $'```postgres-sql\n' >> $REPORT_NAME
cat q3-9.sql >> $REPORT_NAME
cat q10-16.sql >> $REPORT_NAME
echo $'```\n\n' >> $REPORT_NAME

echo $'# PostGIS\n' >> $REPORT_NAME
echo $'```postgres-sql\n' >> $REPORT_NAME
cat q17-19-postgis.sql >> $REPORT_NAME
cat q20-23.sql >> $REPORT_NAME

echo $'```\n' >> $REPORT_NAME
echo "done"
