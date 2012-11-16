import os
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

db_file = os.path.join(os.path.dirname(__file__),'../../_db/seqtool.db')

# echo=True for debug
engine = create_engine('sqlite:///'+db_file, echo=False)
Session = sessionmaker(bind=engine)

