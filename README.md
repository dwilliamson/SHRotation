
##### Spherical Harmonic Rotation
##### Don Williamson (http://donw.io), July 2004
---

I've had this implementation on my HD for quite a while, trying to find time to prepare it for the website. There were a few bugs in the original implementation which Richard Furse was kind enough to point out and I've repackaged everything so as to be useful in realtime game code. This is only a second design pass on the so you'll no doubt want to refactor the code yourself if you'll ever use it.

The simple goal is, provided with a 3x3 rotation matrix, generate an SH rotation matrix, non-analytically, over an abitrary number of bands. Two prominent methods for achieving this are provided by Ivanic/Ruedenberg and Choi. Implementations of both are provided that adhere to the respective papers (and corrections).

NOTE:
Choi's method works by calculating the rotation matrices using complex numbers. We require real rotation matrices and the same paper documents a conversion process from complex to real matrices. Unfortunately the paper contains errors and there is no available corrections from the author. I've spoken to Choi about the errors and he's provided possible solutions (documented in the source code) but a more thourough analysis is required. If anybody has the time to do so, I would love to see the results!


##### Use
---

Original OpenGL rotated spherical function is on the left.
SH rotated spherical function is on the right.

Hold [Space] = match rotation with Ivanic/Ruedenberg  
Hold [Space] + [c] = match rotation with Choi  



##### Acknowledgements
---

Thanks to Richard Furse for highlighing errors in
the original code, Robin Green for the NRook pointer, Cheol Ho Choi for paying
heed to my emails on his original paper, and Joseph Ivanic and Klaus Ruedenberg
for discussing their original algorithm.


##### References
---

In the source code.
