stopAfter=overlapper

# original asm settings
utgErrorRate = 0.25
utgErrorLimit = 4.5
 
cnsErrorRate = 0.25
cgwErrorRate = 0.25
ovlErrorRate = 0.25
 
#merSize=14 # this was here by default
 
# Serge Koren's corrections for 75 bp:
merSize=10
frgMinLen = 30
ovlMinLen = 30
 
 
merylMemory = 128000
merylThreads = 12
 
ovlStoreMemory = 32768
 
# grid info
useGrid = 0
scriptOnGrid = 0
frgCorrOnGrid = 1
ovlCorrOnGrid = 1
 
sge = -A assembly
sgeScript = -pe threads 12
sgeConsensus = -pe threads 1
sgeOverlap = -pe threads 2
sgeFragmentCorrection = -pe threads 2
sgeOverlapCorrection = -pe threads 1
 
#ovlMemory=1GB --hashload 0.7
ovlHashBits=23
ovlHashBlockLength= 30000000
ovlThreads = 2
ovlHashBlockLength = 20000000
ovlRefBlockSize =  5000000
 
# for mer overlapper
merCompression = 1
merOverlapperSeedBatchSize = 500000
merOverlapperExtendBatchSize = 250000
 
frgCorrThreads = 2
frgCorrBatchSize = 100000
 
ovlCorrBatchSize = 100000
 
# non-Grid settings, if you set useGrid to 0 above these will be used
merylMemory = 128000
merylThreads = 12
 
ovlStoreMemory = 32768
 
ovlConcurrency = 50
 
cnsConcurrency = 50
 
merOverlapperThreads = 6
merOverlapperSeedConcurrency = 45
merOverlapperExtendConcurrency = 45
 
frgCorrConcurrency = 50
ovlCorrConcurrency = 50
cnsConcurrency = 50

