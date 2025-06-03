"""
`buffer` has the following structure:

* buffer[1:3,:] contains target positions
* buffer[4:6,:] contains source positions
* buffer[7,:] contains source strength
"""
function direct_buffer!(buffer)
	#--- potential ---#
	
	# 1/r kernel
	@views buffer[1:3,:] .-= buffer[4:6,:]
	@views buffer[4:6,:] .= buffer[1:3]
	@views buffer[4:6,:] .*= buffer[4:6,:]
	@views buffer[4,:] .+= buffer[5,:]
	@views buffer[4,:] .+= buffer[6,:] # r^2
	buffer[5,:] .= 1.0
	@views buffer[5,:] ./= buffer[4,:] # 1/r^2
	@views buffer[4,:] .= sqrt.(buffer[5,:]) # 1/r

	# scale by strength
	@views buffer[4,:] .*= 0.07957747154594767 # 1/r/4pi
	@views buffer[4,:] .*= buffer[7,:] # m/r/4pi (POTENTIAL RESULT)

	#--- velocity ---#

	@views buffer[1:3] .*= buffer[5,:] # dx / r^2
	@views buffer[1:3] .*= buffer[6,:] # dx * m / r^3 / 4pi (VELOCITY RESULT)

end

function fill_direct!(buffer, target_system, target_index, source_system, source_index)
    # get buffer size
    buffer_size = size(buffer, 2)

    
end