function [DSI_gradients] = generate_DSI_vectors(DSI_gradients_number)

% generate DSI vectors from the gradient numbers
% lwkuo, Jan, 2005
% usage: [DSI_gradients] = generate_DSI_vectors(DSI_gradients_number)

DSI_gradients = [];
for p = -10:1:10
    for q = -10:1:10
        for r = -10:1:10
            DSI_gradients = [DSI_gradients;[p^2+q^2+r^2,p,q,r]];
        end
    end
end
DSI_gradients = sortrows(DSI_gradients,1);
DSI_gradients = DSI_gradients(1:DSI_gradients_number,2:4);
