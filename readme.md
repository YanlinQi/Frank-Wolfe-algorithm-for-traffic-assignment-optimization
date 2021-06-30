In this project, the Frank-Wolf algorithm was implemented in Python to solve the traffic assignment problem for the Anaheim network.

* The link performance function for each link is given by the BPR link performance function $t_a = t_{a}^{0}\left [ 1+\alpha \left ( x_a/c_{a}^{'}\right ) ^\beta \right ],
\alpha =0.15, \beta =4, c_{a}^{'}=0.9c_a$.

* The Anaheim network contains 19 origins, 19 destinations (361 O-D pairs), 914 links, and 416 nodes.  

* The topology of the network is given in file anaheim.xlsx.  
* This spreadsheet also contains the travel cost (calculated by free flow speed and the link length) and capacity information of each link. 

* The OD trip table is given in fort2.txt along with a spreadsheet file which explains its format (this is a forward-star representation of the OD information).
