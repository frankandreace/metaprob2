################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BloomFilter/GenomicBloomFilter.cpp \
../src/BloomFilter/RepeatBloomFilter.cpp 

OBJS += \
./src/BloomFilter/GenomicBloomFilter.o \
./src/BloomFilter/RepeatBloomFilter.o 

CPP_DEPS += \
./src/BloomFilter/GenomicBloomFilter.d \
./src/BloomFilter/RepeatBloomFilter.d 


# Each subdirectory must supply rules for building sources it contributes
src/BloomFilter/%.o: ../src/BloomFilter/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


