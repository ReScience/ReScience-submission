"""
   Copyright 2016 Aaron R. Shifman

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

import sys
import model

if sys.argv[1] == "f1":
    model.figure1()
elif sys.argv[1] == "f2":
    model.figure2()
elif sys.argv[1] == "f3":
    model.figure3()
elif sys.argv[1] == "f4":
    model.figure4()
elif sys.argv[1] == "f5":
    model.figure5()
elif sys.argv[1] == "f6":
    model.figure6()
elif sys.argv[1] == "f7":
    model.figure7()
elif sys.argv[1] == "f8":
    model.figure8()
elif sys.argv[1] == "f9":
    model.figure9()
else:
    model.figure1()
    model.figure2()
    model.figure3()
    model.figure4()
    model.figure5()
    model.figure6()
    model.figure7()
    model.figure8()
    model.figure9()