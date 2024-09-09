import sys
from pathlib import Path
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))

from commscoresutil import CommScoresUtil
from scores.antismash import antiSMASH
from scores.bss import bss
from scores.cip import cip
from scores.fs import fs
from scores.gyd import gyd
from scores.mip import mip
from scores.mro import mro
from scores.pc import pc
from scores.smetana import mp, mu, sc, smetana
