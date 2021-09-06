################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/notUsed/Assessment.cpp \
../src/notUsed/CodingSequence.cpp 

OBJS += \
./src/notUsed/Assessment.o \
./src/notUsed/CodingSequence.o 

CPP_DEPS += \
./src/notUsed/Assessment.d \
./src/notUsed/CodingSequence.d 


# Each subdirectory must supply rules for building sources it contributes
src/notUsed/%.o: ../src/notUsed/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


