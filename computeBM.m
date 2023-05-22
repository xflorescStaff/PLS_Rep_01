function UAV_BM = computeBM(dAB, hj, B_vec, rJ_vec )
    
    [~, m_ind] = min(abs(B_vec - dAB));

    UAV_BM = [rJ_vec(m_ind,:), hj];
    
%     fprintf('dAB = %.4f \t\t Bmin = %.4f\n', dAB, B_vec(m_ind))
    
end