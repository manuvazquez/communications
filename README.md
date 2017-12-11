# communications
Source code for different communications-related simulations. Among other things, it encompasses the code for reproducing the results in:

* [_User Activity Tracking in DS-CDMA_][DSCMA], Manuel A. Vázquez, Joaquín Míguez. IEEE Transactions on Vehicular Technology. 62/7, pp. 3188-3203. (*CDMASystem* class)
* [_A Per-Survivor Processing Receiver for
MIMO Transmission Systems With One Unknown Channel Order Per Output_][PSP], Manuel A. Vázquez, Joaquín Míguez. IEEE Transactions on Vehicular Technology. 60/9, pp. 4415-4426. (*ISWCS10System* class)
* [_Maximum-Likelihood Sequence Detection
in Time- and Frequency-Selective MIMO Channels with Unknown Order_][unknown order], Manuel A. Vázquez, Joaquín Míguez. IEEE Transactions on Vehicular
Technology. 58/1, pp. 499-504. (*TVT2007System* and *Rev2TVT2007System* classes)
* [_Sequential Monte Carlo methods for complexity-constrained MAP equalization of dispersive MIMO channels_][complexity constrained], Manuel A. Vázquez, Mónica Fernández Bugallo, Joaquín Míguez. Signal processing. 88/4, pp. 1017-1034. (*Elsevier2007ARChannelSystem* and *Elsevier2007BesselChannelSystem* classes)
* [_On the Use of the Channel Second-Order
Statistics in MMSE receivers for time- and frequency-selective MIMO transmission systems_][SOS], Manuel A. Vázquez, Joaquín Míguez. EURASIP Journal on
Wireless Communications and Networking. (*PlainSystem* class)


[DSCMA]: http://ieeexplore.ieee.org/abstract/document/6473922/
[PSP]: http://ieeexplore.ieee.org/document/6032763/
[unknown order]: http://ieeexplore.ieee.org/abstract/document/4510724/
[complexity constrained]: http://www.sciencedirect.com/science/article/pii/S0165168407003763
[SOS]: https://jwcn-eurasipjournals.springeropen.com/articles/10.1186/s13638-016-0768-0

Requirements
============

- [eigen](http://eigen.tuxfamily.org/)
- [rapidxml](http://rapidxml.sourceforge.net/)

The code has only been tested when compiled with gcc.
