#!/bin/sh

if [ ! -n "$4" ]
then
  echo "Usage: `basename $0` dir host user[:password] db"
  exit 65
fi  

DIR=$1
HOST=$2
USERPASS=$3
DB=$4

USER=${USERPASS%:*}
echo ${USER}
PASS=${USERPASS#*:}
echo ${PASS}

cd $DIR
for SQL in *.sql
do
    T_NAME=${SQL%%.sql}
    echo "loading table ${T_NAME}"
    if [ $USERPASS == $USER ]
    then
	mysql --compress --host=$HOST --user=$USER -e "DROP TABLE IF EXISTS ${T_NAME};" ${DB} > /dev/null 2> /dev/null
	mysql --compress --host=$HOST --user=$USER ${DB} < ${SQL}
	zcat "${T_NAME}.txt.gz" | mysql --compress --host=$HOST --user=$USER --local-infile=1 \
	    -e "LOAD DATA LOCAL INFILE \"/dev/stdin\" INTO TABLE ${T_NAME};" ${DB}
    else
	mysql --compress --host=$HOST --user=$USER --password=$PASS -e "DROP TABLE IF EXISTS ${T_NAME};" ${DB} > /dev/null 2> /dev/null
	mysql --compress --host=$HOST --user=$USER --password=$PASS ${DB} < ${SQL}
	zcat "${T_NAME}.txt.gz" | mysql --compress --host=$HOST --user=$USER --password=$PASS --local-infile=1 \
	    -e "LOAD DATA LOCAL INFILE \"/dev/stdin\" INTO TABLE ${T_NAME};" ${DB}
    fi
done
