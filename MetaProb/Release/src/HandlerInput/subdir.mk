################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/HandlerInput/DNAFile.cpp \
../src/HandlerInput/Input.cpp \
../src/HandlerInput/Parameter.cpp 

OBJS += \
./src/HandlerInput/DNAFile.o \
./src/HandlerInput/Input.o \
./src/HandlerInput/Parameter.o 

CPP_DEPS += \
./src/HandlerInput/DNAFile.d \
./src/HandlerInput/Input.d \
./src/HandlerInput/Parameter.d 


# Each subdirectory must supply rules for building sources it contributes
src/HandlerInput/%.o: ../src/HandlerInput/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


