classdef TryClass
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        meshrefiner
        easyvisualizer
        udmwriter
        nodepos
        shearrates
        temps
        fillstatus
        connectednodes
        connectedelements
        velx
        vely
        velz
    end

    methods
        function obj = ObjectConstuctor(obj, orig_mesh)
            addpath('..\..\utilities');
            obj.meshrefiner = meshRefiner;
            obj.easyvisualizer = easyVisualizer;
            obj.udmwriter = udmWriter;
            obj.nodepos = load(fullfile('input', orig_mesh, 'Np.mat')).Np;
            obj.shearrates = load(fullfile('input', orig_mesh, 'Shearrate.mat')).SRnode;
            obj.shearrates(obj.shearrates(:, :) == -1e+30) = 0;
            obj.temps = load(fullfile('input', orig_mesh, 'Tnode.mat')).Tnode;
            obj.fillstatus = load(fullfile('input', orig_mesh, 'Fillstatus.mat')).Fillstatus;
            obj.connectednodes = load(fullfile('input', orig_mesh, 'connectedNodes.mat')).connectedNodes;
            obj.connectedelements = load(fullfile('input', orig_mesh, 'connectedElements.mat')).connectedElements;
            obj.velx = load(fullfile('input', orig_mesh, 'VelX.mat')).VelX;
            obj.vely = load(fullfile('input', orig_mesh, 'VelY.mat')).VelY;
            obj.velz = load(fullfile('input', orig_mesh, 'VelZ.mat')).VelZ;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end