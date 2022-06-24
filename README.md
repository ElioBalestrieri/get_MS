# get_MS
MS detection algorithm


implementation of MS detection algorithm described in Engbert & Kliegl (2003) with further modifications in Engbert and Mergenthaler *PNAS* (2006).

Function optimized to work on data coming from Buschlab even if the function is created for handling data from 2 eyes, 
it does not contain an algorithm of selection of MS computed from both eyes simultaneously.

It is advisable to average the signal from both eyes before submitting it to the function (as in Lin, Nobre and Van Ede, Nature Comms 2022).

