Forecasting the incidence of jamming attacks in Wireless Sensor Networks using Epidemic Logistic Growth Model
Author: Miguel Lopez, PhD
BioSiP Research Group - University of Malaga - Spain 
e-mail. m.lopez@uma.es // code.io@icloud.com
Revision date.  Dec. 12, 2022


Scenario 1: Reached nodes by the coordinator when the jammer node is near to the network's coordinator
Scenario 2: Reached nodes by the coordinator when the jammer node is in the middle of the topology
Scenario 3: Reached nodes by the coordinator when the jammer node is far to the network's coordinator
Type of jamming: Random & Reactive
Protocols involved: IEEE 802.15.4, AODV, DSR & MPH
Cases of study from 1 to 3 are 50 p/s random jamming
               from 4 to 6 are 80 p/s random jamming
               from 7 to 9 are reactive jamming attack

NOTES:

1. Place in the same folder that the MATLAB code the jamming data files named rawjammingdata.txt

2. Run the GGM model ggm_jamming to obtain a first set of estimated values for r and p, at the early phase of the jamming attack

3. Run the GLGM model glgm_jamming to obtain the forecasted curves of the attack, and the values for r, p, K and Ro of the jamming attack.

