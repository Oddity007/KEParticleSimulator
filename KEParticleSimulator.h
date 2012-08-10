/*
Copyright (c) 2012 Oliver Daids

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

	1. The origin of this software must not be misrepresented; you must not
	claim that you wrote the original software. If you use this software
	in a product, an acknowledgment in the product documentation would be
	appreciated but is not required.

	2. Altered source versions must be plainly marked as such, and must not be
	misrepresented as being the original software.

	3. This notice may not be removed or altered from any source
	distribution.
*/

#ifndef KEParticleSimulator_h
#define KEParticleSimulator_h
#ifdef __cplusplus
extern "C"
{
#endif
#include <stdint.h>

/*typedef uint8_t KEParticleSimulatorTypeCode;
#define KEParticleSimulatorTypeCodeNone 0
#define KEParticleSimulatorTypeCodeUInt8 1
#define KEParticleSimulatorTypeCodeUInt16 2
#define KEParticleSimulatorTypeCodeUInt32 3
#define KEParticleSimulatorTypeCodeUInt64 4
#define KEParticleSimulatorTypeCodeInt8 5
#define KEParticleSimulatorTypeCodeInt16 6
#define KEParticleSimulatorTypeCodeInt32 7
#define KEParticleSimulatorTypeCodeInt64 8
#define KEParticleSimulatorTypeCodeFloat32 11
#define KEParticleSimulatorTypeCodeFloat64 12*/

typedef struct KEParticleSimulator KEParticleSimulator;
typedef uint32_t KEParticleSimulatorClusterID;
typedef enum
{
	KEParticleSimulatorClusterTypeRadialForce,
	KEParticleSimulatorClusterTypeDirectionalForce,
	KEParticleSimulatorClusterTypeParticle,
	
}KEParticleSimulatorClusterType;
/*typedef enum
{
	KEParticleSimulatorClusterPropertyNameCenters,
	KEParticleSimulatorClusterPropertyNameRadii,
	KEParticleSimulatorClusterPropertyNameVelocities,
	KEParticleSimulatorClusterPropertyNameMasses,
	KEParticleSimulatorClusterPropertyNameRadialForces,
	KEParticleSimulatorClusterPropertyNameDirectionalForces,
	
	KEParticleSimulatorClusterPropertyNameInteractionMask
}KEParticleSimulatorClusterPropertyName;*/

KEParticleSimulator* KEParticleSimulatorNew(void);
void KEParticleSimulatorDelete(KEParticleSimulator* self);

const char* KEParticleSimulatorPollErrorMessages(KEParticleSimulator* self);

void KEParticleSimulatorUpdate(KEParticleSimulator* self, double seconds);

KEParticleSimulatorClusterID KEParticleSimulatorCreateCluster(KEParticleSimulator* self, KEParticleSimulatorClusterType type);
void KEParticleSimulatorDestroyCluster(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID);

void KEParticleSimulatorSetClusterBounds(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID, const double center[3], const double halfBounds[3]);

/*void KEParticleSimulatorBeginClusterEditing(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID);
void KEParticleSimulatorSetClusterPropertyPointer(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID, KEParticleSimulatorClusterPropertyName propertyName, uint32_t elementStart, uint32_t elementCount, uint32_t elementComponentCount, KEParticleSimulatorTypeCode elementComponentTypeCode, bool isNormalized, size_t stride, const void* pointer);
void KEParticleSimulatorEndClusterEditing(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID);

void KEParticleSimulatorGetClusterPropertyPointer(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID, KEParticleSimulatorClusterPropertyName propertyName, uint32_t elementStart, uint32_t elementCount, uint32_t elementComponentCount, KEParticleSimulatorTypeCode elementComponentTypeCode, bool isNormalized, size_t stride, void* pointer);*/

//This method keeps things simple
typedef struct KEParticleSimulatorClusterElementRadialForce KEParticleSimulatorClusterElementRadialForce;
typedef struct KEParticleSimulatorClusterElementDirectionalForce KEParticleSimulatorClusterElementDirectionalForce;
typedef struct KEParticleSimulatorClusterElementParticle KEParticleSimulatorClusterElementParticle;

struct KEParticleSimulatorClusterElementRadialForce
{
	float
		position[3],
		force,
		radius;
};

struct KEParticleSimulatorClusterElementDirectionalForce
{
	float
		position[3],
		force[3],
		radius;
};

struct KEParticleSimulatorClusterElementParticle
{
	float
		position[3],
		velocity[3],
		mass;
};

typedef uint64_t KEParticleSimulatorClusterElementMappingMode;
enum
{
	KEParticleSimulatorClusterElementMappingModeNone = 0,
	KEParticleSimulatorClusterElementMappingModeRead = (1 << 0),
	KEParticleSimulatorClusterElementMappingModeWrite = (1 << 1),
};

void* KEParticleSimulatorMapClusterElements(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID, uint32_t elementCount, KEParticleSimulatorClusterElementMappingMode mappingMode);
void KEParticleSimulatorUnmapClusterElements(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID);

#ifdef __cplusplus
}
//extern "C"
#endif

#endif
