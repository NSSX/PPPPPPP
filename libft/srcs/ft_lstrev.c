/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_lstrev.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/24 19:15:15 by igomez            #+#    #+#             */
/*   Updated: 2014/11/27 12:19:06 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

void	ft_lstrev(t_list **lst_ptr)
{
	t_list	*prev;
	t_list	*temp;

	if (!lst_ptr || !(*lst_ptr))
		return ;
	prev = NULL;
	temp = NULL;
	while (temp || !prev)
	{
		temp = (*lst_ptr)->next;
		(*lst_ptr)->next = prev;
		prev = *lst_ptr;
		if (temp)
			*lst_ptr = temp;
	}
}
