extern "C"
{
#include "KEParticleSimulator.h"
}
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>

struct Cluster
{
	KEParticleSimulatorClusterType type;
	void* elements;
	uint32_t elementCount;
	bool isAllocated;
	KEParticleSimulatorClusterElementMappingMode mappingMode;
	double overdueTimeLeft;
};

struct KEParticleSimulator
{
	std::vector<Cluster> clusters;
};

KEParticleSimulator* KEParticleSimulatorNew(void)
{
	KEParticleSimulator* self = new KEParticleSimulator;
	return self;
}

void KEParticleSimulatorDelete(KEParticleSimulator* self)
{
	if(not self) return;
	for (KEParticleSimulatorClusterID i = 0; i < self->clusters.size(); i++)
	{
		Cluster* cluster = & self->clusters[i];
		if(not cluster->isAllocated) continue;
		free(cluster->elements);
	}
	delete self;
}

const char* KEParticleSimulatorPollErrorMessages(KEParticleSimulator* self)
{
	return NULL;
}


void KEParticleSimulatorUpdate(KEParticleSimulator* self, double seconds)
{
	const float secondsPerUpdate = 1.0f / 30.0f;//Eventually make this a user changeable value
	for (KEParticleSimulatorClusterID i = 0; i < self->clusters.size(); i++)
	{
		Cluster* cluster = & self->clusters[i];
		if(not cluster->isAllocated) continue;
		switch (cluster->type)
		{
			case KEParticleSimulatorClusterTypeParticle:
				cluster->overdueTimeLeft += seconds;
				while(cluster->overdueTimeLeft > secondsPerUpdate)
				{
					for (uint32_t j = 0; j < cluster->elementCount; j++)
					{
						KEParticleSimulatorClusterElementParticle* particle = j + (KEParticleSimulatorClusterElementParticle*) cluster->elements;
						//Begin brute-force, unoptimized force accumulation
						for (KEParticleSimulatorClusterID i = 0; i < self->clusters.size(); i++)
						{
							Cluster* forceCluster = & self->clusters[i];
							if(not forceCluster->isAllocated) continue;
							switch (forceCluster->type)
							{
								case KEParticleSimulatorClusterTypeParticle:
									break;
								case KEParticleSimulatorClusterTypeDirectionalForce:
									for (uint32_t k = 0; k < forceCluster->elementCount; k++)
									{
										KEParticleSimulatorClusterElementDirectionalForce* directionalForce = k + (KEParticleSimulatorClusterElementDirectionalForce*) forceCluster->elements;
										float distanceSquared = 0;
										for (int l = 0; l < 3; l++)
										{
											float difference = particle->position[l] - directionalForce->position[l];
											distanceSquared += difference * difference;
										}
										if(distanceSquared > (directionalForce->radius * directionalForce->radius)) continue;
										float inverseDistance = 1.0f/sqrtf(distanceSquared);
										//We multiply by the square inverse so we can normalize and attenuate at the same time
										float accelerationFactor = (1 / particle->mass) * inverseDistance * inverseDistance * secondsPerUpdate;
										for (int l = 0; l < 3; l++)
											particle->velocity[l] += directionalForce->force[l] * accelerationFactor;
									}
									break;
								case KEParticleSimulatorClusterTypeRadialForce:
									for (uint32_t k = 0; k < forceCluster->elementCount; k++)
									{
										KEParticleSimulatorClusterElementRadialForce* radialForce = k + (KEParticleSimulatorClusterElementRadialForce*) forceCluster->elements;
										float distanceSquared = 0;
										float direction[3];
										for (int l = 0; l < 3; l++)
										{
											float difference = particle->position[l] - radialForce->position[l];
											direction[l] = difference;
											distanceSquared += difference * difference;
										}
										if(distanceSquared > (radialForce->radius * radialForce->radius)) continue;
										float inverseDistance = 1.0f/sqrtf(distanceSquared);
										//We multiply by the square inverse so we can normalize and attenuate at the same time
										float acceleration = (radialForce->force / particle->mass) * inverseDistance * inverseDistance * secondsPerUpdate;
										for (int l = 0; l < 3; l++)
											particle->velocity[l] += direction[l] * acceleration;
									}
									break;
							}
						}
						//End force accumulation
						for (int k = 0; k < 3; k++)
							particle->position[k] += particle->velocity[k] * secondsPerUpdate;
					}
					cluster->overdueTimeLeft -= secondsPerUpdate;
				}
				break;
			default:
				break;
		}
	}
}

KEParticleSimulatorClusterID KEParticleSimulatorCreateCluster(KEParticleSimulator* self, KEParticleSimulatorClusterType type, uint32_t elementCount)
{
	Cluster cluster;
	cluster.isAllocated = true;
	cluster.type = type;
	cluster.elementCount = elementCount;
	switch (type)
	{
		case KEParticleSimulatorClusterTypeDirectionalForce:
			cluster.elements = calloc(elementCount, sizeof(KEParticleSimulatorClusterElementDirectionalForce));
			break;
		case KEParticleSimulatorClusterTypeParticle:
			cluster.elements = calloc(elementCount, sizeof(KEParticleSimulatorClusterElementParticle));
			break;
		case KEParticleSimulatorClusterTypeRadialForce:
			cluster.elements = calloc(elementCount, sizeof(KEParticleSimulatorClusterElementRadialForce));
			break;
	}
	
	for (KEParticleSimulatorClusterID i = 0; i < self->clusters.size(); i++)
	{
		if(self->clusters[i].isAllocated) continue;
		self->clusters[i] = cluster;
		return i + 1;
	}
	
	self->clusters.push_back(cluster);
	return self->clusters.size();
}

void KEParticleSimulatorDestroyCluster(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID)
{
	if(not clusterID) return;
	Cluster* cluster = & self->clusters[clusterID - 1];
	free(cluster->elements);
	cluster->isAllocated = false;
	if(self->clusters.size() == clusterID) self->clusters.pop_back();
}

void* KEParticleSimulatorMapClusterElements(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID, KEParticleSimulatorClusterElementMappingMode mappingMode)
{
	if(not clusterID) return NULL;
	Cluster* cluster = & self->clusters[clusterID - 1];
	cluster->mappingMode = mappingMode;
	return cluster->elements;
}

void KEParticleSimulatorUnmapClusterElements(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID)
{
	if(not clusterID) return;
	Cluster* cluster = & self->clusters[clusterID - 1];
	KEParticleSimulatorClusterElementMappingMode mappingMode = cluster->mappingMode;
	cluster->mappingMode = KEParticleSimulatorClusterElementMappingModeNone;
	if(not (mappingMode & KEParticleSimulatorClusterElementMappingModeWrite)) return;
	//Update cache or whatever
}
