rsync -aP --delete database/ ubuntu@54.93.53.183:~/database
ssh ubuntu@54.93.53.183 'sudo supervisorctl restart all'

