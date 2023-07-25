# install packages
pip3 install -r ./requirements.txt

# run API server
## localhost: gunicorn mhw_app:app -k uvicorn.workers.UvicornWorker -b 127.0.0.1:8020
gunicorn mhw_app:app -w 4 -k uvicorn.workers.UvicornWorker -b 127.0.0.1:8020 --keyfile conf/privkey.pem --certfile conf/fullchain.pem --reload

# kill process
ps -ef | grep 'gunicorn' | grep -v grep | awk '{print $2}' | xargs -r kill -9

# pm2 start
pm2 start ./conf/ecosystem.config.js

