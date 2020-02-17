                function [ spikeMat , tVec ] = poissonSpikeGen ( fr , tSim , nTrials, dt )
                
                %{
                http://www.columbia.edu/cu/appliedneuroshp/Fall2017/neuralcoding.pdf
                
                %}
                
                %dt = 0.001; % s
                nBins = floor ( tSim / dt ) ;
                spikeMat = rand ( nTrials , nBins ) < fr * dt ;
                tVec = 0: dt : tSim - dt ;
                end