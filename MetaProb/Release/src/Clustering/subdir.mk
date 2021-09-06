################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Clustering/ClusterUtility.cpp \
../src/Clustering/Gmeans.cpp \
../src/Clustering/KStest.cpp \
../src/Clustering/Kmeans.cpp 

OBJS += \
./src/Clustering/ClusterUtility.o \
./src/Clustering/Gmeans.o \
./src/Clustering/KStest.o \
./src/Clustering/Kmeans.o 

CPP_DEPS += \
./src/Clustering/ClusterUtility.d \
./src/Clustering/Gmeans.d \
./src/Clustering/KStest.d \
./src/Clustering/Kmeans.d 


# Each subdirectory must supply rules for building sources it contributes
src/Clustering/%.o: ../src/Clustering/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


