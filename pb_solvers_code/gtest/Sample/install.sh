#!/bin/bash

rm *.a
make clean
make gtest-all.o
make gtest_main.o
make gtest.a
make gtest_main.a
make testsuite
