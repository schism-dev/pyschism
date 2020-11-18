import os
from typing import List

from flask import Flask
from flask_apscheduler import APScheduler  # type: ignore[import]
from waitress import serve

# from flask_restful import Resource, Api
# from sqlalchemy import create_engine  # type: ignore[import]

#  def get_sqlite_engine(echo=False):
#     return create_engine(f'sqlite:///{str(FORECASTD_DATA_DIR)}/index.sqlite',
#                          echo=echo)


class Config:
    JOBS: List = []
    SCHEDULER_API_ENABLED = True


# class JobList(Resource):

#     def get(self, name):
#         return {'job_name': name}

app = Flask(__name__)
app.config.from_object(Config())
scheduler = APScheduler()
scheduler.init_app(app)
scheduler.start()
# api = Api(app)
# scheduler = BackgroundScheduler(get_sqlite_engine())
# print()
# self._scheduler.add_job(run_forecast_cycle, 'interval', seconds=30)
# self._scheduler.start()


@app.route('/')
def hello():
    return 'HELLO!!!'
    # for i in range(10):
    #     app.apscheduler.add_job(func=scheduled_task, trigger='date', args=[i], id='j'+str(i))

# return 'Scheduled several long running tasks.', 200


def run():
    serve(
        app,
        host=os.getenv(
            'FORECASTD_HOSTNAME', 'http://localhost').split('://')[-1],
        port=os.getenv('FORECASTD_PORT', 5000),
        # debug=True
    )
