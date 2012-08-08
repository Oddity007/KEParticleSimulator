#include "KEParticleSimulator.h"
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>

namespace
{
	struct Cluster
	{
		KEParticleSimulatorClusterType type;
		void* elements;
		uint32_t elementCount;
		bool isAllocated;
		KEParticleSimulatorClusterElementMappingMode mappingMode;
		double overdueTimeLeft;
	};
}

struct KEParticleSimulator
{
	std::vector<Cluster> clusters;
};

extern "C"
{

KEParticleSimulator* KEParticleSimulatorNew(void)
{
	KEParticleSimulator* self = new KEParticleSimulator;
	return self;
}

void KEParticleSimulatorDelete(KEParticleSimulator* self)
{
	if(not self) return;
	//Clean up any remaining allocations
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
	//No error messages implemented yet
	return NULL;
}

//This update loop's code is temporary until I performance test it and optimize.
void KEParticleSimulatorUpdate(KEParticleSimulator* self, double seconds)
{
	//Eventually make this a user changeable value
	const float secondsPerUpdate = 1.0f / 30.0f;
	//For each cluster
	for (std::vector<Cluster>::iterator cluster = self->clusters.begin(); cluster not_eq self->clusters.end(); cluster++)
	//for (uint32_t i = 0; i < self->clusters.size(); i++)
	{
		//Cluster* cluster = &(self->clusters[i]);
		//Check if it is in use or not, skip to the next one if it isn't
		if(not cluster->isAllocated) continue;
		switch (cluster->type)
		{
			case KEParticleSimulatorClusterTypeParticle:
				//Lock to a fixed update rate of secondsPerUpdate
				cluster->overdueTimeLeft += seconds;
				while(cluster->overdueTimeLeft >= secondsPerUpdate)
				{
					for (uint32_t j = 0; j < cluster->elementCount; j++)
					{
						KEParticleSimulatorClusterElementParticle* particle = j + (KEParticleSimulatorClusterElementParticle*) cluster->elements;
						//Begin brute-force, unoptimized force accumulation
						//Iterate over each cluster again
						for (std::vector<Cluster>::iterator forceCluster = self->clusters.begin(); forceCluster not_eq self->clusters.end(); forceCluster++)
						{
							if(not forceCluster->isAllocated) continue;
							switch (forceCluster->type)
							{
								case KEParticleSimulatorClusterTypeParticle:
									//Skip particle clusters
									break;
								case KEParticleSimulatorClusterTypeDirectionalForce:
									for (uint32_t k = 0; k < forceCluster->elementCount; k++)
									{
										const KEParticleSimulatorClusterElementDirectionalForce* directionalForce = k + (KEParticleSimulatorClusterElementDirectionalForce*) forceCluster->elements;
										//Calcuate the distance squared between the particle and force
										float distanceSquared = 0;
										for (int l = 0; l < 3; l++)
										{
											float difference = particle->position[l] - directionalForce->position[l];
											distanceSquared += difference * difference;
										}
										//Skip this force if the particle is out of range
										if(distanceSquared > (directionalForce->radius * directionalForce->radius)) continue;
										//We need the inverse distance for attenuation
										float inverseDistance = 1.0f / sqrtf(distanceSquared);
										//Attenuate by the inverse distance, scale the force by the inverse mass, and scale by the update time
										float accelerationFactor = (1.0f / particle->mass) * inverseDistance * secondsPerUpdate;
										//Now add to the velocity
										for (int l = 0; l < 3; l++)
											particle->velocity[l] += directionalForce->force[l] * accelerationFactor;
									}
									break;
								case KEParticleSimulatorClusterTypeRadialForce:
									for (uint32_t k = 0; k < forceCluster->elementCount; k++)
									{
										const KEParticleSimulatorClusterElementRadialForce* radialForce = k + (KEParticleSimulatorClusterElementRadialForce*) forceCluster->elements;
										//Calcuate the distance squared and direction vector between the particle and force.  We don't normalize the direction just yet though.
										float distanceSquared = 0;
										float direction[3];
										for (int l = 0; l < 3; l++)
										{
											float difference = particle->position[l] - radialForce->position[l];
											direction[l] = difference;
											distanceSquared += difference * difference;
										}
										//Skip this force if the particle is out of range
										if(distanceSquared > (radialForce->radius * radialForce->radius)) continue;
										//We need the inverse distance for attenuation
										float inverseDistance = 1.0f/sqrtf(distanceSquared);
										//Attenuate by the inverse distance, scale the force by the inverse mass, and scale by the update time
										//We multiply by the square inverse distance so we can normalize and attenuate at the same time
										float acceleration = (radialForce->force / particle->mass) * inverseDistance * inverseDistance * secondsPerUpdate;
										//Now add to the velocity
										for (int l = 0; l < 3; l++)
											particle->velocity[l] += direction[l] * acceleration;
									}
									break;
							}
						}
						//End force accumulation
						//Now update the particle position with the velocity
						for (int k = 0; k < 3; k++)
							particle->position[k] += particle->velocity[k] * secondsPerUpdate;
					}
					//Decrease the timer
					cluster->overdueTimeLeft -= secondsPerUpdate;
					//The process will repeat again until cluster->overdueTimeLeft < secondsPerUpdate
				}
				break;
			default:
				break;
		}
	}
}

KEParticleSimulatorClusterID KEParticleSimulatorCreateCluster(KEParticleSimulator* self, KEParticleSimulatorClusterType type)
{
	Cluster cluster;
	cluster.isAllocated = true;
	cluster.type = type;
	cluster.elementCount = 0;
	cluster.elements = NULL;
	
	//Find an empty slot if it exists
	for (KEParticleSimulatorClusterID i = 0; i < self->clusters.size(); i++)
	{
		if(self->clusters[i].isAllocated) continue;
		self->clusters[i] = cluster;
		return i + 1;
	}
	
	//Otherwise, make a new one
	self->clusters.push_back(cluster);
	return (KEParticleSimulatorClusterID) self->clusters.size();
}

void KEParticleSimulatorDestroyCluster(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID)
{
	if(not clusterID) return;
	Cluster* cluster = & self->clusters[clusterID - 1];
	free(cluster->elements);
	cluster->isAllocated = false;
	if(self->clusters.size() == clusterID) self->clusters.pop_back();
}

void* KEParticleSimulatorMapClusterElements(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID, uint32_t elementCount, KEParticleSimulatorClusterElementMappingMode mappingMode)
{
	if(not clusterID) return NULL;
	Cluster* cluster = & self->clusters[clusterID - 1];
	cluster->mappingMode = mappingMode;
	//if it's the same size, we just exit
	if(cluster->elementCount == elementCount) return cluster->elements;
	//otherwise, we need to do a reallocation
	size_t elementSize = 0;
	switch (cluster->type)
	{
		case KEParticleSimulatorClusterTypeDirectionalForce:
			elementSize = sizeof(KEParticleSimulatorClusterElementDirectionalForce);
			break;
		case KEParticleSimulatorClusterTypeParticle:
			elementSize = sizeof(KEParticleSimulatorClusterElementParticle);
			break;
		case KEParticleSimulatorClusterTypeRadialForce:
			elementSize = sizeof(KEParticleSimulatorClusterElementRadialForce);
			break;
	}
	cluster->elements = realloc(cluster->elements, elementSize * elementCount);
	//If the new allocation is larger than the old allocation, zero out the new memory regions, but leave the old data in the old regions
	if(cluster->elementCount < elementCount) memset(((char*)cluster->elements) + elementSize * cluster->elementCount, 0, elementSize * (elementCount - cluster->elementCount));
	//And finally update the element count to reflect the new one
	cluster->elementCount = elementCount;
	return cluster->elements;
}

void KEParticleSimulatorUnmapClusterElements(KEParticleSimulator* self, KEParticleSimulatorClusterID clusterID)
{
	if(not clusterID) return;
	Cluster* cluster = & self->clusters[clusterID - 1];
	KEParticleSimulatorClusterElementMappingMode mappingMode = cluster->mappingMode;
	cluster->mappingMode = KEParticleSimulatorClusterElementMappingModeNone;
	if(not (mappingMode & KEParticleSimulatorClusterElementMappingModeWrite)) return;
	//Update cache, or whatever, here
}

}
//extern  "C"