#!/bin/bash

./main ./data/DeLucia2006a_80M 2 ./out_files/DeLucia2006a_80M_2 &
./main ./data/DeLucia2006a_80M 4 ./out_files/DeLucia2006a_80M_4 &
./main ./data/DeLucia2006a_80M 8 ./out_files/DeLucia2006a_80M_8 &
./main ./data/DeLucia2006a_80M 16 ./out_files/DeLucia2006a_80M_16 &
./main ./data/DeLucia2006a_80M 32 ./out_files/DeLucia2006a_80M_32 &
./main ./data/DeLucia2006a_80M 64 ./out_files/DeLucia2006a_80M_64 &
wait
./main ./data/DeLucia2006a_100M 2 ./out_files/DeLucia2006a_100M_2 &
./main ./data/DeLucia2006a_100M 4 ./out_files/DeLucia2006a_100M_4 &
./main ./data/DeLucia2006a_100M 8 ./out_files/DeLucia2006a_100M_8 &
./main ./data/DeLucia2006a_100M 16 ./out_files/DeLucia2006a_100M_16 &
./main ./data/DeLucia2006a_100M 32 ./out_files/DeLucia2006a_100M_32 &
./main ./data/DeLucia2006a_100M 64 ./out_files/DeLucia2006a_100M_64 &
wait