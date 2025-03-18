#!/bin/bash
# Run on local machine to connect to kopernik

ssh -J aszary@kepler.ia.uz.zgora.pl aszary@kopernik.ia.uz.zgora.pl -o ServerAliveInterval=60 -o ServerAliveCountMax=5

