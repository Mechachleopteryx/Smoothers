## -*- encoding: utf-8 -*-
# This file was *autogenerated* from the file smoothers.sagetex.sage.
from sage.all_cmdline import *   # import sage library
_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_1p5 = RealNumber('1.5'); _sage_const_28 = Integer(28); _sage_const_16 = Integer(16); _sage_const_20 = Integer(20); _sage_const_26 = Integer(26)## This file (smoothers.sagetex.sage) was *autogenerated* from smoothers.tex with sagetex.sty version 2012/01/16 v2.3.3-69dcb0eb93de.
import sagetex
_st_ = sagetex.SageTeXProcessor('smoothers', version='2012/01/16 v2.3.3-69dcb0eb93de', version_check=True)
_st_.blockbegin()
try:
 sage.misc.preparser.load(sage.misc.preparser.base64.b64decode("YmFzaWMuc2FnZQ=="),globals(),False)
 patterns=[[-_sage_const_3 ,_sage_const_0 ],[_sage_const_1 ,-_sage_const_0 ],[_sage_const_1p5 ,-_sage_const_0 ]]
 g1=shl()
except:
 _st_.goboom(_sage_const_16 )
_st_.blockend()
try:
 _st_.plot(_sage_const_0 , format='notprovided', _p_=g1)
except:
 _st_.goboom(_sage_const_20 )
_st_.blockbegin()
try:
 f=lfunctions(_sage_const_1 /_sage_const_2 )[_sage_const_1 ].simplify()
except:
 _st_.goboom(_sage_const_26 )
_st_.blockend()
try:
 _st_.inline(_sage_const_0 , latex(f))
except:
 _st_.goboom(_sage_const_28 )
_st_.endofdoc()
