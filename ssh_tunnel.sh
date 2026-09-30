#!/bin/bash
# Run on local machine to enable VNC connection

# Check if a port number is provided
if [ -z "$1" ]; then
	echo "Usage: $0 <port-number e.g 5900>"
	exit 1
fi

PORT=$1

ssh -L $PORT:kopernik:$PORT -J aszary@kepler.ia.uz.zgora.pl aszary@kopernik.ia.uz.zgora.pl -o ServerAliveInterval=60 -o ServerAliveCountMax=5

