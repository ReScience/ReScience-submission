'''
This script was written for the purpose of extracting the plotted points in Fig. 2 and Fig. 4 of Riehle et al. (1997)
Science 278, 1950-1953.
The extraction should be performed as follows:
1) prepare a copy of the paper PDF.
2) generate an EPS file containing only the figure part of the PDF file as follows:
    2-1) open the PDF file in Adobe Acrobat Pro (the authors used Adobe Acrobat XI Pro on Windows 7 OS).
    2-2) toggle text editing by clicking "Edit">"Edit Text & Images".
    2-3) delete everything in the document except for the figure of concern.
    2-4) export the result to an EPS file.
3) Run this script after setting the variable `filename` to the name of the generated EPS file (without extension).
The extracted points are saved in a NumPY .npy format file with the same name as the EPS file.
Note: the positions of the extracted points are represented in arbitrary units defined in the original PDF file.

Copyright (c) 2016-2017, Vahid Rostami, Junji ito, Michael Denker, Sonja Gruen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
following conditions are met:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following
disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
following disclaimer in the documentation and/or other materials provided with the distribution.
Neither the names of the copyright holders nor the names of the contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
import numpy as np


filename = "Riehle et al._1997_Science_fig2"
color_dict = {
    (1, 0, 0, 0): "cyan",
    (0, 1, 0, 0): "magenta",
    (0, 0, 1, 0): "yellow",
    (1, 1, 0, 0): "blue",
    (0, 1, 1, 0): "red",
    (1, 0, 1, 0): "green",
}
event_type_dict = {
    "black": 0,  # spike
    "cyan": 1,  # coincidence
    "red": 2,  # UE
}

xs_tmp = []
ys_tmp = []
xs = []
ys = []
event_types = []
mo_found = False
color = "black"

# A brief explanation of the logic:
# The main loop below scans lines in the given EPS file and looks for a pair of successive lines looking like:
#
#     *x* *y* mo
#     *x* *y* li
#
# where *x* and *y* are arbitrary numbers.
# These lines constitute a command set for drawing a point at position (*x*, *y*).
# At every occurrence of this pattern, *x* and *y* are extracted and kept in temporal buffer lists `xs_tmp` and
# `ys_tmp`.
# The loop also looks for a line specifying the color of drawing, i.e., a line ending with "cmyk", and read the color
# information (referreing to the `color_dict` dictionary defined above) to store it in a temporal buffer variable
# `color`.
# The collected position and color information is copied to outcome buffer lists `xs`, `ys` and `event_types`, once
# a line composed only of a single "@", which is a command to execute the drawing command sets requested so far, is
# encountered during the scan.
# Note that the color information is converted to event type information before stored in the list, referring to the
# `event_type_dict` dictionary defined above.
# After storing the information in the outcome buffers, the temporal buffer lists are flushed and the scanning of
# lines starts over.

for line in open(filename + ".eps", 'r'):
    token = line.split()
    if len(token) == 1:
        if token[0] == "@" and len(xs_tmp) > 0:
            xs.extend(xs_tmp)
            ys.extend(ys_tmp)
            event_types.extend([event_type_dict[color]] * len(xs_tmp))
            xs_tmp = []
            ys_tmp = []
        mo_found = False
    elif len(token) == 3:
        if token[-1] == "mo":
            x_pre, y_pre = float(token[0]), float(token[1])
            mo_found = True
        elif token[-1] == "li":
            x, y = float(token[0]), float(token[1])
            if mo_found and x == x_pre and y == y_pre:
                xs_tmp.append(x)
                ys_tmp.append(-y)
            mo_found = False
        else:
            mo_found = False
    elif len(token) == 5:
        if token[-1] == "cmyk":
            c, m, y, k = [int(x) for x in token[:4]]
            color = color_dict[(c, m, y, k)]
        else:
            color = "black"
        mo_found = False
    elif len(token) > 5:
        color = "black"
        mo_found = False
    else:
        mo_found = False

data = np.array(zip(xs, ys, event_types),
                dtype=[('x', np.float), ('y', np.float),
                       ('event_type', np.int)])
np.save(filename + ".npy", data)

# an example usage of the saved data
data_saved = np.load(filename + ".npy")
import matplotlib.pyplot as plt
plt.plot(data_saved['x'][data_saved['event_type'] == 0],
         data_saved['y'][data_saved['event_type'] == 0], 'k.')
plt.plot(data_saved['x'][data_saved['event_type'] == 1],
         data_saved['y'][data_saved['event_type'] == 1], 'c.')
plt.plot(data_saved['x'][data_saved['event_type'] == 2],
         data_saved['y'][data_saved['event_type'] == 2], 'm.')
plt.show()
