function GC3D_plot_seg_results(run_param, seg_param)
        
   %% Plotting the result of segmentation
        figure(run_param.fignum);
        isosurface(run_param.Label1);
        %         h = gcf;
        %         dump_PDF_figure(h, [seg_param.sample_name '_' seg_param.algo_type '_dist']);

        
    %% Function to get draw nice contours on the object 
        function out = get_contours(cc)
        s = 1;
        k = 2;
        out = {};
            while k <= size(cc, 2)
                m = cc(2, k - 1);
                out{s} = cc(:, k:(k+m-1));
                k = k + m + 1;
                s = s + 1;
            end 