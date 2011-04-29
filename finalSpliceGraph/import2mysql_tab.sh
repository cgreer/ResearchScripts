#!/bin/bash

if [ ! -n "$5" ]
then
  echo "Usage: `basename $0` dir table host user[:password] db"
  exit 65
fi  

DIR=$1
T_NAME=$2
HOST=$3
USERPASS=$4
DB=$5

USER=${USERPASS%:*}
PASS=${USERPASS#*:}

cd $DIR
if [ -e "$T_NAME.sql" ]
then
    echo "loading table ${T_NAME}"
    if [ $USERPASS == $USER ]
    then
	mysql --compress --host=$HOST --user=$USER -e "DROP TABLE IF EXISTS ${T_NAME};" ${DB} > /dev/null 2> /dev/null
	echo "here one"
	mysql --compress --host=$HOST --user=$USER ${DB} < ${T_NAME}.sql
	echo "here"
	zcat "${T_NAME}.txt.gz" | mysql --compress --host=$HOST --user=$USER --local-infile=1 \
	    -e "LOAD DATA LOCAL INFILE \"/dev/stdin\" INTO TABLE ${T_NAME};" ${DB}
	echo "finished importing data into table ${T_NAME}"
    else
	echo "second route: ${PASS}"
	mysql --compress --host=$HOST --user=$USER --password=$PASS -e "DROP TABLE IF EXISTS ${T_NAME};" ${DB} > /dev/null 2> /dev/null
	mysql --compress --host=$HOST --user=$USER --password=$PASS ${DB} < ${T_NAME}.sql
	zcat "${T_NAME}.txt.gz" | mysql --compress --host=$HOST --user=$USER --password=$PASS --local-infile=1 \
	    -e "LOAD DATA LOCAL INFILE \"/dev/stdin\" INTO TABLE ${T_NAME};" ${DB}
    fi
else
    echo "skipped table: $T_NAME.sql does not exist in $DIR"
fi
