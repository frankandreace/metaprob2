################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Cluster.cpp \
../src/CountNucleotide.cpp \
../src/Group.cpp \
../src/HashTable.cpp \
../src/Normalization.cpp \
../src/SeqGraph.cpp \
../src/SequencePE.cpp \
../src/Solution.cpp \
../src/Utilities.cpp \
../src/VirSeq.cpp \
../src/main.cpp 

OBJS += \
./src/Cluster.o \
./src/CountNucleotide.o \
./src/Group.o \
./src/HashTable.o \
./src/Normalization.o \
./src/SeqGraph.o \
./src/SequencePE.o \
./src/Solution.o \
./src/Utilities.o \
./src/VirSeq.o \
./src/main.o 

CPP_DEPS += \
./src/Cluster.d \
./src/CountNucleotide.d \
./src/Group.d \
./src/HashTable.d \
./src/Normalization.d \
./src/SeqGraph.d \
./src/SequencePE.d \
./src/Solution.d \
./src/Utilities.d \
./src/VirSeq.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


