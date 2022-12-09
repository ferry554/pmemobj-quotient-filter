# pmemobj-quotient-filter

A NVM-based quotient filter written in C,using libpmemobj.

pmem-qf.c: Implementation  
pmem-qf.h: API and documentation  
test.cc: Randomized tester  

To build:  
`make test`  

To verify its correctness:  
`./test <filename> test`  

To measure the performance of insert and lookup:  
`./test <filename> bench`  
