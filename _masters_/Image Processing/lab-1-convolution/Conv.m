classdef Conv
    properties (GetAccess = private)
        kernel_x;
        kernel_y;
    end
    
    properties
        name;
    end
    
    methods
        function obj = Conv(m_name, m_x, m_y)
            obj.name = m_name;
            obj.kernel_x = m_x;
            obj.kernel_y = m_y;
        end
        
        function convolute(obj, img)
            w_x = conv2(img, obj.kernel_x); 
            w_y = conv2(img, obj.kernel_y);
            w = sqrt(double(w_x.^2 + w_y.^2)); % wiki formula
            
            % draw image
            filename = [obj.name '.png'];
            imwrite(w, gray(256), filename);

            % draw convolution plot
            figure();
            mesh(w);
            filename = [obj.name '-plot.png'];
            saveas(gcf, filename);
        end
    end
end
