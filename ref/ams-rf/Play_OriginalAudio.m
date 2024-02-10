[InputAudio,fsAudio]    = audioread('lab_data/out.wav');
InputAudio              = InputAudio(1:2^19,1)/max(InputAudio(:,1));
soundsc(InputAudio,fsAudio);