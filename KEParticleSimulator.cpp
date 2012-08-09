#include "KEParticleSimulator.h"
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>

#pragma mark Utility

namespace
{
	static inline float InvSqrt(float x)
	{
		return 1.0f / sqrtf(x);
		/*float xhalf = 0.5f * x;
		int i = *(int*)&x; // store floating-point bits in integer
		i = 0x5f3759d5 - (i >> 1); // initial guess for Newton's method
		x = *(float*)&i; // convert new bits into float
		x = x*(1.5f - xhalf*x*x); // One round of Newton's method
		return x;*/
	}

	struct AABB
	{
		float center[3], halfBounds[3];
		
		//Checks if "this" touches the supplied parameters
		static bool CheckIntersection(const AABB& a, const AABB& b)
		{
			for (int i = 0; i < 3; i++)
			{
				float difference = fabsf(a.center[i] - b.center[i]);
				float boundSum = a.halfBounds[i] + b.halfBounds[i];
				if(difference >= boundSum) return false;
			}
			return true;
		}
		
		//Checks if a encapsulates the supplied b
		static bool CheckEncapsulation(const AABB& a, const AABB& b)
		{
			for (int i = 0; i < 3; i++)
			{
				float aMax = a.center[i] + a.halfBounds[i];
				float bMax = b.center[i] + b.halfBounds[i];
				if(aMax < bMax) return false;
				float aMin = a.center[i] - a.halfBounds[i];
				float bMin = b.center[i] - b.halfBounds[i];
				if(aMin > bMin) return false;
			}
			return true;
		}
		
		static AABB Union(const AABB& a, const AABB& b)
		{
			AABB result;
			for (int i = 0; i < 3; i++)
			{
				float aMax = a.center[i] + a.halfBounds[i];
				float bMax = b.center[i] + b.halfBounds[i];
				float aMin = a.center[i] - a.halfBounds[i];
				float bMin = b.center[i] - b.halfBounds[i];
				float max = (aMax > bMax) ? aMax : bMax;
				float min = (aMin > bMin) ? aMin : bMin;
				result.halfBounds[i] = (max - min) * 0.5f;
				result.center[i] = min + result.halfBounds[i];
			}
			return result;
		}
	};
}

#pragma mark Particle Simulation

namespace
{
	struct Cluster
	{
		KEParticleSimulatorClusterType type;
		void* elements;
		uint32_t elementCount;
		bool isAllocated;
		KEParticleSimulatorClusterElementMappingMode mappingMode;
		AABB bounds;
	};
	
	struct UpdateableClusterCache
	{
		std::vector<const Cluster*> intersectingDirectionalForceClusters;
		std::vector<const Cluster*> intersectingRadialForceClusters;
		struct Package
		{
			Cluster* particleCluster;
			uint32_t
				intersectingDirectionalForceClusterCount,
				intersectingRadialForceClusterCount;
		};
		std::vector<Package> updateableParticleClusterPackages;
	};
}

struct KEParticleSimulator
{
	std::vector<Cluster> clusters;
	UpdateableClusterCache* updateableClusterCache;
	double overdueTimeLeft;
};

namespace
{
	static void RecalculateSingleSimulatorClusterCache(KEParticleSimulator* self, Cluster* cluster)
	{
		switch (cluster->type)
		{
			case KEParticleSimulatorClusterTypeParticle:
				{
					const KEParticleSimulatorClusterElementParticle* elements = (const KEParticleSimulatorClusterElementParticle*) cluster->elements;
					float
						minimums[3],
						maximums[3];
					memcpy(minimums, elements[0].position, sizeof(minimums));
					memcpy(maximums, elements[0].position, sizeof(maximums));
					for (uint32_t i = 0; i < cluster->elementCount; i++)
					{
						for (int j = 0; j < 3; j++)
						{
							float position = elements[i].position[j];
							maximums[j] = (position > maximums[j]) ? position : maximums[j];
							minimums[j] = (position < minimums[j]) ? position : minimums[j];
						}
					}
					for (int i = 0; i < 3; i++)
					{
						cluster->bounds.halfBounds[i] = (maximums[i] - minimums[i]) * 0.5f;
						cluster->bounds.center[i] = cluster->bounds.halfBounds[i] + minimums[i];
					}
				}
				break;
			case KEParticleSimulatorClusterTypeDirectionalForce:
				{
					const KEParticleSimulatorClusterElementDirectionalForce* elements = (const KEParticleSimulatorClusterElementDirectionalForce*) cluster->elements;
					float
						minimums[3],
						maximums[3];
					memcpy(minimums, elements[0].position, sizeof(minimums));
					memcpy(maximums, elements[0].position, sizeof(maximums));
					for (uint32_t i = 0; i < cluster->elementCount; i++)
					{
						for (int j = 0; j < 3; j++)
						{
							float maximum = elements[i].position[j] + elements[i].radius;
							maximums[j] = (maximum > maximums[j]) ? maximum : maximums[j];
							float minimum = elements[i].position[j] - elements[i].radius;
							minimums[j] = (minimum < minimums[j]) ? minimum : minimums[j];
						}
					}
					for (int i = 0; i < 3; i++)
					{
						cluster->bounds.halfBounds[i] = (maximums[i] - minimums[i]) * 0.5f;
						cluster->bounds.center[i] = cluster->bounds.halfBounds[i] + minimums[i];
					}
				}
				break;
			case KEParticleSimulatorClusterTypeRadialForce:
				{
					const KEParticleSimulatorClusterElementRadialForce* elements = (const KEParticleSimulatorClusterElementRadialForce*) cluster->elements;
					float
						minimums[3],
						maximums[3];
					memcpy(minimums, elements[0].position, sizeof(minimums));
					memcpy(maximums, elements[0].position, sizeof(maximums));
					for (uint32_t i = 0; i < cluster->elementCount; i++)
					{
						for (int j = 0; j < 3; j++)
						{
							float maximum = elements[i].position[j] + elements[i].radius;
							maximums[j] = (maximum > maximums[j]) ? maximum : maximums[j];
							float minimum = elements[i].position[j] - elements[i].radius;
							minimums[j] = (minimum < minimums[j]) ? minimum : minimums[j];
						}
					}
					for (int i = 0; i < 3; i++)
					{
						cluster->bounds.halfBounds[i] = (maximums[i] - minimums[i]) * 0.5f;
						cluster->bounds.center[i] = cluster->bounds.halfBounds[i] + minimums[i];
					}
				}
				break;
		}
	}
	
	static void RegenerateSimulatorUpdateableClusterCache(KEParticleSimulator* self)
	{
		//if(self->updateableClusterCache) return;
		if(not self->updateableClusterCache) self->updateableClusterCache = new UpdateableClusterCache;
		//Clear out the old data
		self->updateableClusterCache->intersectingDirectionalForceClusters.clear();
		self->updateableClusterCache->intersectingRadialForceClusters.clear();
		self->updateableClusterCache->updateableParticleClusterPackages.clear();
		//For each cluster
		for (uint32_t particleClusterIndex = 0; particleClusterIndex < self->clusters.size(); particleClusterIndex++)
		//for (uint32_t i = 0; i < self->clusters.size(); i++)
		{
			Cluster* particleCluster = &(self->clusters[particleClusterIndex]);
			UpdateableClusterCache::Package package;
			//Check if it is in use or not, skip to the next one if it isn't
			if(not particleCluster->isAllocated) continue;
			//Skip if it's not a particle cluster
			if(particleCluster->type not_eq KEParticleSimulatorClusterTypeParticle) continue;
			//If we have no elements, skip
			if(not particleCluster->elementCount) continue;
			//Now collect all force clusters that intersect this cluster
			package.intersectingDirectionalForceClusterCount = 0;
			package.intersectingRadialForceClusterCount = 0;
			package.particleCluster = particleCluster;
			for (uint32_t forceClusterIndex = 0; forceClusterIndex < self->clusters.size(); forceClusterIndex++)
			{
				const Cluster* forceCluster = &(self->clusters[forceClusterIndex]);
				if(not forceCluster->isAllocated) continue;
				if(not AABB::CheckIntersection(forceCluster->bounds, particleCluster->bounds)) continue;
				switch (forceCluster->type)
				{
					case KEParticleSimulatorClusterTypeParticle: break;
					case KEParticleSimulatorClusterTypeDirectionalForce:
						self->updateableClusterCache->intersectingDirectionalForceClusters.push_back(forceCluster);
						package.intersectingDirectionalForceClusterCount++;
						break;
					case KEParticleSimulatorClusterTypeRadialForce:
						self->updateableClusterCache->intersectingRadialForceClusters.push_back(forceCluster);
						package.intersectingRadialForceClusterCount++;
						break;
				}
			}
			self->updateableClusterCache->updateableParticleClusterPackages.push_back(package);
		}
	}

	static void ClearSimulatorUpdateableClusterCache(KEParticleSimulator* self)
	{
		delete self->updateableClusterCache;
		self->updateableClusterCache = NULL;
	}
	
	static void UpdateParticleVelocityWithDirectionalForce(KEParticleSimulatorClusterElementParticle* particle, const KEParticleSimulatorClusterElementDirectionalForce* directionalForce, double secondsPerUpdate)
	{
		//Calcuate the distance squared between the particle and force
		float distanceSquared = 0;
		for (int l = 0; l < 3; l++)
		{
			float difference = particle->position[l] - directionalForce->position[l];
			distanceSquared += difference * difference;
		}
		//Skip this force if the particle is out of range
		if(distanceSquared > (directionalForce->radius * directionalForce->radius)) return;
		//We need the inverse distance for attenuation
		float inverseDistance = InvSqrt(distanceSquared);
		//Attenuate by the inverse distance, scale the force by the inverse mass, and scale by the update time
		float accelerationFactor = (1.0f / particle->mass) * inverseDistance * secondsPerUpdate;
		//Now add to the velocity
		for (int l = 0; l < 3; l++)
			particle->velocity[l] += directionalForce->force[l] * accelerationFactor;
	}

	static void UpdateParticleVelocityWithRadialForce(KEParticleSimulatorClusterElementParticle* particle, const KEParticleSimulatorClusterElementRadialForce* radialForce, double secondsPerUpdate)
	{
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
		if(distanceSquared > (radialForce->radius * radialForce->radius)) return;
		//We need the inverse distance for attenuation
		float inverseDistance = InvSqrt(distanceSquared);
		//Attenuate by the inverse distance, scale the force by the inverse mass, and scale by the update time
		//We multiply by the square inverse distance so we can normalize and attenuate at the same time
		float acceleration = (radialForce->force / particle->mass) * inverseDistance * inverseDistance * secondsPerUpdate;
		//Now add to the velocity
		for (int l = 0; l < 3; l++)
			particle->velocity[l] += direction[l] * acceleration;
	}

	static void UpdateSimulator(KEParticleSimulator* self, double secondsPerUpdate)
	{
		uint32_t
			intersectingDirectionalForceClusterIterationBaseIndex = 0,
			intersectingRadialForceClusterIterationBaseIndex = 0;
	
		//Iterate the packages
		for (uint32_t packageIndex = 0; packageIndex < self->updateableClusterCache->updateableParticleClusterPackages.size(); packageIndex++)
		{
			const UpdateableClusterCache::Package& package = self->updateableClusterCache->updateableParticleClusterPackages[packageIndex];
			{
				for (uint32_t j = 0; j < package.particleCluster->elementCount; j++)
				{
					KEParticleSimulatorClusterElementParticle* particle = j + (KEParticleSimulatorClusterElementParticle*) package.particleCluster->elements;
					//Begin brute-force, unoptimized force accumulation
					//Iterate over each intersecting directional force cluster
					for (uint32_t i = 0; i < package.intersectingDirectionalForceClusterCount; i++)
					{
						const Cluster* forceCluster = self->updateableClusterCache->intersectingDirectionalForceClusters[intersectingDirectionalForceClusterIterationBaseIndex + i];
						for (uint32_t k = 0; k < forceCluster->elementCount; k++)
						{
							const KEParticleSimulatorClusterElementDirectionalForce* directionalForce = k + (KEParticleSimulatorClusterElementDirectionalForce*) forceCluster->elements;
							UpdateParticleVelocityWithDirectionalForce(particle, directionalForce, secondsPerUpdate);
						}
					}
				
					//Iterate over each intersecting radial force cluster
					for (uint32_t i = 0; i < package.intersectingRadialForceClusterCount; i++)
					{
						const Cluster* forceCluster = self->updateableClusterCache->intersectingRadialForceClusters[intersectingRadialForceClusterIterationBaseIndex + i];
						for (uint32_t k = 0; k < forceCluster->elementCount; k++)
						{
							const KEParticleSimulatorClusterElementRadialForce* radialForce = k + (KEParticleSimulatorClusterElementRadialForce*) forceCluster->elements;
							UpdateParticleVelocityWithRadialForce(particle, radialForce, secondsPerUpdate);
						}
					}
					//End force accumulation
					//Now update the particle position with the velocity
					for (int k = 0; k < 3; k++)
							particle->position[k] += particle->velocity[k] * secondsPerUpdate;
				}
			}
			//Shift the intersecting cluster iteration base indices
			intersectingDirectionalForceClusterIterationBaseIndex += package.intersectingDirectionalForceClusterCount;
			intersectingRadialForceClusterIterationBaseIndex += package.intersectingRadialForceClusterCount;
		}
		
		//Now regenerate the AABBs.  This is done separately, as we can get away with doing this less often if we want.
		for (uint32_t packageIndex = 0; packageIndex < self->updateableClusterCache->updateableParticleClusterPackages.size(); packageIndex++)
		{
			const UpdateableClusterCache::Package& package = self->updateableClusterCache->updateableParticleClusterPackages[packageIndex];
			RecalculateSingleSimulatorClusterCache(self, package.particleCluster);
		}
	}
}

extern "C"
{

KEParticleSimulator* KEParticleSimulatorNew(void)
{
	KEParticleSimulator* self = new KEParticleSimulator;
	self->updateableClusterCache = NULL;
	self->overdueTimeLeft = 0;
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
	ClearSimulatorUpdateableClusterCache(self);
	delete self;
}

const char* KEParticleSimulatorPollErrorMessages(KEParticleSimulator* self)
{
	//No error messages implemented yet
	return NULL;
}

void KEParticleSimulatorUpdate(KEParticleSimulator* self, double seconds)
{
	RegenerateSimulatorUpdateableClusterCache(self);

	//Eventually make this a user changeable value
	const float secondsPerUpdate = 1.0f / 30.0f;
	
	//Lock to a fixed update rate of secondsPerUpdate
	self->overdueTimeLeft += seconds;
	while(self->overdueTimeLeft >= secondsPerUpdate)
	{
		UpdateSimulator(self, seconds);
		//Decrease the timer
		self->overdueTimeLeft -= secondsPerUpdate;
		//The process will repeat again until self->overdueTimeLeft < secondsPerUpdate
	}
}

KEParticleSimulatorClusterID KEParticleSimulatorCreateCluster(KEParticleSimulator* self, KEParticleSimulatorClusterType type)
{
	ClearSimulatorUpdateableClusterCache(self);

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
	ClearSimulatorUpdateableClusterCache(self);
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
	ClearSimulatorUpdateableClusterCache(self);
	RecalculateSingleSimulatorClusterCache(self, cluster);
}

}
//extern  "C"