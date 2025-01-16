import csv
import os

from Data.Analyze.tools.batch.compile_logs import get_logs_and_pdbs
from Data.Analyze.tools.compare.read_logs2 import read_logs2
import pandas as pd
from System.sys_funcs.calcs.calcs import calc_dist
from System.system import System

